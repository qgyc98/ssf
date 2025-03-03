#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "galias.h"
#include "iotools.h"
#include "prepalias.h"
#include "stacktrace.h"

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
long Edg;

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
  xfscanf(in, "%a", Fname);
  //xfscanf(in, "%a", Afname);
  Nn = (Nexo+1) * (Neyo+1);
  Ne = Nexo * Neyo*2;
  fprintf(stdout, "\nNumber of nodes : %ld", Nn);
  fprintf(stdout, "\nNumber of elems : %ld", Ne);
  fprintf(stdout, "\nFile name for output : %s\n", Fname);
  //xfscanf(in, "%a", Tmpl);
  //Edg = 0;
  //fscanf(in, "%ld", &Edg);
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
  long did;
  El = new long *[Ne];
  memset(El, 0, sizeof(*El)*Ne);
  Len = new long *[Ndx*Ndy];
  memset(Len, 0, sizeof(*Len)*Ndx*Ndy*2);

  for (i = 0; i < Ndx; i++)
  {
    for (j = 0; j < Ndy; j++)
    {
      lie = 0;
      did = i * Ndy + j;
      Len[did] = new long[2*Nex[i]*Ney[j]];
      memset(Len[did], 0, sizeof(*Len[did])*Nex[i]*Ney[j]);
      for (ii = 0; ii < Nex[i]; ii++){
        for (jj = 0; jj < Ney[j]; jj++){
          El[ie] = new long[3];
          memset(El[ie], 0, sizeof(*El[ie])*3);
	  // first element in rectangle
          El[ie][0] = (ii+1)*(Ney[j]+1)+jj+1;
	  El[ie][1] = ii*(Ney[j]+1)+jj;
	  El[ie][2] = (ii+1)*(Ney[j]+1)+jj;
	  Len[did][lie] = ie;
          lie++;
          ie++;
	  //	  the second element in rectangle
	  El[ie] = new long[3];
	  memset(El[ie], 0, sizeof(*El[ie])*3);
	  El[ie][0] = (ii+1)*(Ney[j]+1)+jj+1;
	  El[ie][1] = ii*(Ney[j]+1)+jj+1;
	  El[ie][2] = ii*(Ney[j]+1)+jj;
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
  long *edge_prop;
  
  edge_prop = new long[3];
  //memset(edge_prop, 0, sizeof(*edge_prop)*3);
  
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
	    //TKo
	    if ((i == 0) && (ii == 0)) // first column of domains and first column of nodes
	    {
	      if ((j == 0) && (jj == 0))
            {
              // bottom left corner has vertex property 1, edge property 1 and 4, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 1  2 1  2 4  3 %ld  4 1\n", did+1);
              id++;
              continue;
            }
            if ((j == Ndy-1) && (jj == Ney[j]))
            {
              // top left corner has vertex property 4, edge property 4 and 3, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 4  2 3  2 4  3 %ld  4 1\n", did+1);
              id++;
              continue;
            }
            // nodes on the left edge have edge property 4, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 4  3 %ld  4 1\n", did+1);
            id++;
            continue;
          }
          if ((i == Ndx-1) && (ii == Nex[i])) // last column of domains and last column of nodes
          {
            if ((j == 0) && (jj == 0))
            {
              // bottom right corner has vertex property 2, edge property 1 and 2, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 2  2 1  2 2  3 %ld  4 1\n", did+1);
              id++;
              continue;
            }
            if ((j == Ndy-1) && (jj == Ney[j]))
            {
              // top right corner has vertex property 3, edge property 2 and 3, surface property 
              // by domain id and volume property 1
              fprintf(out, "5  1 3  2 2  2 3  3 %ld  4 1\n", did+1);
              id++;
              continue;
            }
            // nodes on the right edge have edge property 2, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 2  3 %ld  4 1\n", did+1);
            id++;
            continue;
          }
          if ((j == 0) && (jj == 0))
          {
            // internal nodes on the bottom edge have edge property 1, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 1  3 %ld  4 1\n", did+1);
            id++;
            continue;
          }
          if ((j == Ndy-1) && (jj == Ney[j]))
          {
            // internal nodes on the top edge have edge property 3, surface property by 
            // domain id and volume property 1
            fprintf(out, "3  2 3  3 %ld  4 1\n", did+1);
            id++;
            continue;
          }
          // internal nodes have surface property by domain id and volume property 1 
          fprintf(out, "2  3 %ld  4 1\n", did+1);
          id++;
	  }
      }

      //JB
      // PRINTING OF ELEMENTS
      id = 0;
      fprintf(out, "%ld\n", Nex[i]*Ney[j]*2);
      for (ii = 0; ii < Nex[i]; ii++)
      {
        for (jj = 0; jj < 2*Ney[j]; jj++)
	  {
          fprintf(out, "\n%ld 3", id+1);
          for (k = 0; k < 3; k++){
            fprintf(out, " %ld", El[Len[did][id]][k]+1);
	    id++;
	  }
          fprintf(out, " 1"); // element volume property
	  
          if (Edg)
          {
	    // type of domain 0
	    if(i == 0 && j == 0){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if(ii >= 0 && jj == 0){
		edge_prop[0] = 0;
		edge_prop[1] = 1;
		edge_prop[2] = 0;
	      }
	      if(ii == 0 && (jj > 0 && ((jj+1)%2 == 0))){
		edge_prop[0] = 0;
		edge_prop[1] = 4;
		edge_prop[2] = 0;
	      } 
	      id++;
	    }
	    //type of domain 1
	    if(i == 0 && (j > 0 && j < Ndy - 1) ){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if(ii == 0 && (jj >= 0 && (jj+1)%2 == 0)){
		edge_prop[0] = 0;
		edge_prop[1] = 4;
		edge_prop[2] = 0;
	      } 
	      id++;
	    }
	    
	    // type of domain 2
	    if(i == 0 && j == Ndy-1){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if(ii == 0 && ((jj > 0 && jj < 2*Ney[j]-1) && (jj+1)%2 == 0 )){
		edge_prop[0] = 0;
		edge_prop[1] = 4;
		edge_prop[2] = 0;
	      }
	    if(ii == 0 && jj == 2*Ney[j]-1){
		edge_prop[0] = 3;
		edge_prop[1] = 4;
		edge_prop[2] = 0;
	      }
	    if(ii > 0 && jj == 2*Ney[j]-1){
		edge_prop[0] = 3;
		edge_prop[1] = 0;
		edge_prop[2] = 0;
	      }
	    id++;
	    }
	      // type of domain 3
	    if((i > 0 && i < Ndx-1) && j == 0){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if(ii >= 0 && jj == 0){
		edge_prop[0] = 0;
		edge_prop[1] = 1;
		edge_prop[2] = 0;
	      }
	      id++;
	    } 
	    // type of domain 4 - inner domain
	    if((i > 0 && i < Ndx-1) && (j > 0 && j < Ndy-1)){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      id++;
	    } 
	    // type of domain 5 
	    if((i > 0 && i < Ndx-1) && j == Ndy-1){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if(ii >= 0 && jj == 2*Ney[j]-1){
		edge_prop[0] = 0;
		edge_prop[1] = 3;
		edge_prop[2] = 0;
	      }
	      id++;
	    } 
	    // type of domain 6
	    if(i == Ndx-1 && j == 0){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if((ii >= 0 && ii < Nex[i]-1) && jj == 0){
		edge_prop[0] = 0;
		edge_prop[1] = 1;
		edge_prop[2] = 0;
	      }
	      if(ii == Nex[i]-1 && jj== 0){
		edge_prop[0] = 0;
		edge_prop[1] = 1;
		edge_prop[2] = 2;
	      }
	      if(ii == Nex[i]-1 && (jj > 0 && jj%2 == 0)){
		edge_prop[0] = 0;
		edge_prop[1] = 0;
		edge_prop[2] = 2;
	      }
	      id++;
	    } 
	    
	    
	    
	    // type of domain 7
	    if(i == Ndx-1 && (j > 0 && j < Ndy-1 )){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if(ii == Nex[i]-1 && (jj >= 0 && jj%2 == 0)){
		edge_prop[0] = 0;
		edge_prop[1] = 0;
		edge_prop[2] = 2;
	      }
	      id++;
	    }
	    
	    // type of domain 8
	    if(i == Ndx-1 && j == Ndy-1){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	      if(ii >= 0 && jj == 2*Ney[j]-1){
	      edge_prop[0] = 3;
	      edge_prop[1] = 0;
	      edge_prop[2] = 0;
	    }
	      if(ii == Nex[i]-1 && (jj >= 0 && jj%2 == 0)){
	      edge_prop[0] = 0;
	      edge_prop[1] = 0;
	      edge_prop[2] = 2;
	      }
	      id++;
	    }
	    fprintf(out, "   %ld %ld %ld %ld\n",edge_prop[0],edge_prop[1],edge_prop[2],did+1);     
            
	    
	  }
        }
      }
      //delete []edge_prop;
      
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
      print_err("Section with input files has not been detected", __FILE__, __LINE__, __func__);
      abort();
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
      // reading of input files 
      xf_setsec(tmplf, bsec_str[begsec_files]);
      id_fs = tmplf->id_sec;
      // topology file name
      fprintf(out, "%s\n", bsec_str[begsec_files].alias);
      fprintf(out, "%s\n", topfn);
      
      // material dbase
      // reading of line with material database file name
      xfscanf(tmplf, " %a", str);
      // write line to the preprocesor file
      fprintf(out, "%s\n", str);
      
      // cross section dbase
      // reading of line with cross-section database file name
      xfscanf(tmplf, " %a", str);
      // write line to the preprocesor file
      fprintf(out, "%s\n", str);

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
      //  mesh format (this is sifel format, see GEFEL/galias.h), edges
      fprintf(out, "mesh_format sifel\n");
      fprintf(out, "edge_numbering %ld\n", Edg);

      // end of file section  
      fprintf(out, "%s\n\n", esec_str[begsec_files].alias);

      // probdesc section
      xf_setsec(tmplf, bsec_str[begsec_probdesc]);
      id_pd = tmplf->id_sec;
      fprintf(out, "%s\n", bsec_str[begsec_probdesc].alias);
      fprintf(out, "domain_number %ld   #domain id", did+1);// domain number
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
  //  FILE *out = fopen("vystup","w");
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
  printtop(Fname);
  xfclose(in);
}

