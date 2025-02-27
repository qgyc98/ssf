#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "siftop.h"
#include "iotools.h"


// keywords for beginning of preprocessor sections
enum  bsec {begsec_etype=0, begsec_nodes=1, begsec_elems=2};
const enumstr bsec_str[] = {{"!element type",0}, {"!n",1}, {"!E",2}};

// keywords for ends of preprocessor sections
const kwdset bsec_kwdset(sizeof(bsec_str)/sizeof(*bsec_str), bsec_str);

enum  esec {endsec_etype=0, endsec_nodes=1, endsec_elems=2};
const enumstr esec_str[] = {{"!n",0}, {"!E",1}, {"EPLOT",2}};
const kwdset esec_kwdset(sizeof(esec_str)/sizeof(*esec_str), esec_str);


int main(int argc, char *argv[])
{
  XFILE *in;
  FILE *out;
  siftop *top;
  char *aptr;
  long n_et, nid;
  long *etypes;
  long i, j;
  gtypel gte;


  fprintf(stdout, "\n\n--- Convertor of topology from ANSYS to SIFEL ---\n");
  fprintf(stdout, "-------------------------------------------------\n\n");

  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : ansys2sif input_file_name output_file name\n\n");
    return(1);
  }

  fprintf(stdout, "Opening the ANSYS file %s ...", argv[1]);
  in = xfopen(argv[1],"r");
  if (in == NULL){
    print_err("cannot find ANSYS file %s", __FILE__, __LINE__, __func__, argv[1]);
    return(1);
  }
  in->warning = 1;
  //in->kwdmode = sequent_mode;
  in->kwdmode = sect_mode_fwd;
  in->ignorecase = 0; // case sensitive searching of keywords
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, "Parsing sections of ANSYS file ...");
  // detection of sections in the preprocesor file
  xfdetect_sect(in, bsec_kwdset, esec_kwdset);
  fprintf(stdout, " O.K.\n");

  // reading of element types
  fprintf(stdout, "Reading of element types used in ANSYS ...");
  if (xf_setsec(in, bsec_str[begsec_etype]))
  {
    print_err("cannot find section with element type by keywords '%s', '%s'", 
              __FILE__, __LINE__, __func__, bsec_str[begsec_etype], esec_str[endsec_etype]);
    return(2);
  }
  in->kwdmode = sect_mode_full;
  // detect number of occurences of keyword "ET"
  n_et = getkwd_sect(in, NULL, aptr, "ET", 1);
  if (n_et > 1)
  {
    print_err("number of element types > 1,\n heterogeneous meshes are not supported", 
              __FILE__, __LINE__, __func__);
    return(2);
  }
  etypes = new long[n_et];
  memset(etypes, 0, sizeof(*etypes)*n_et);
  for (i=0; i<n_et; i++)
  {
    getkwd_sect(in, NULL, aptr, "ET", 0);
    xfscanf(in, "%k , %*ld , %ld", "ET", etypes+i);
  }
  switch(etypes[0])
  {
    case 187: 
      gte = tetrahedronlinear;
      break;
    default:
      print_err("unknown type of element is required (et = %d)", 
                __FILE__, __LINE__, __func__, etypes[0]);
      return(2);
  }
  fprintf(stdout, " O.K.\n\n");

  top = new siftop;

  // reading of nodes
  fprintf(stdout, "Importing of nodes ...");
  if (xf_setsec(in, bsec_str[begsec_nodes]))
  {
    print_err("cannot find section with nodes marked by keywords '%s', '%s'", 
              __FILE__, __LINE__, __func__, bsec_str[begsec_nodes], esec_str[endsec_nodes]);
    return(3);
  }
  // detect number of occurences of keyword "n"
  in->kwdmode = sect_mode_full;
  top->nn = getkwd_sect(in, NULL, aptr, "n", 1);
  top->nodes = new snode[top->nn];

  // reading of nodes
  in->warning = 0;
  in->kwdmode = sect_mode_fwd;
  for (i=0; i<top->nn; i++)
  {
    getkwd_sect(in, NULL, aptr, "n", 0);
    xfscanf(in, "%k", "n");
    xfscanf(in, " , %ld", &nid);
    if ((nid < 1) && (nid > top->nn))
    {
      print_err("node number %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)\n",
                __FILE__, __LINE__, __func__, nid, top->nn, in->fname, in->line, in->col);
      return (4);
    }
    nid--;
    xfscanf(in, " , %le , %le , %le", &top->nodes[nid].x, &top->nodes[nid].y, &top->nodes[nid].z);
    top->nodes[nid].alloc(1);
    top->nodes[nid].entid[0] = eregion, 
    top->nodes[nid].prop[0] = 0;
  }
  in->warning = 1;
  fprintf(stdout, "    %8ld of nodes was found\n", top->nn);
  

  // reading of elements
  fprintf(stdout, "Importing of elements ...");
  if (xf_setsec(in, bsec_str[begsec_elems]))
  {
    print_err("cannot find section with nodes marked by keywords '%s', '%s'", 
              __FILE__, __LINE__, __func__, bsec_str[begsec_elems], esec_str[endsec_elems]);
    return(5);
  }
  top->ne = getkwd_sect(in, NULL, aptr, "E", 1);
  top->elements = new selement[top->ne];
  in->warning = 0;
  in->kwdmode = sect_mode_fwd;
  for (i=0; i<top->ne; i++)
  {
    top->elements[i].type = gte;
    if (top->elements[i].alloc(0))
    {
      print_err("unknown type on element %ld is required (file %s, line=%ld, col=%ld)", 
                __FILE__, __LINE__, __func__, i+1, in->fname, in->line, in->col);
      return (6);
    }
    getkwd_sect(in, NULL, aptr, "E", 0);
    xfscanf(in, "%k", "E");
    for (j = 0; j < top->elements[i].nne; j++)
    {
      xfscanf(in, " , %ld", &top->elements[i].nodes[j]);
      top->elements[i].nodes[j]--;
    }
  }
  fprintf(stdout, " %8ld of elements was found\n\n", top->ne);

  in->warning = 1;

  // export of topology in sifel format
  fprintf(stdout, "Opening the SIFEL file %s ...", argv[2]);
  out = fopen(argv[2],"w");
  if (out == NULL){
    fprintf (stderr,"\n Cannot open the output file '%s'.", argv[2]);
    return(1);
  }
  fprintf(stdout, " O.K.\n");
  fprintf(stdout, "Writing topology in SIFEL format ...");
  top->print(out);
  fprintf(stdout, " O.K.\n\n");

  delete top;
  delete [] etypes;
  return(0);
}
