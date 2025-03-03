#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// GEFEL header files
#include "gtopology.h"
#include "siftop.h"
#include "nodmap.h"

// generator header files
#include "lshape_geom.h"
#include "cprof_geom.h"

int main(int argc, char *argv[])
{
  XFILE *in = NULL;
  FILE *out = NULL;
  XFILE *intf = NULL;
  const char *infn;
  const char *outfn;
  infn = "cg-input.in";
  outfn = "cg.top";
  char outhnfn[] = "hangnod.in";
  lshape_geom cg;
  nodmap *cpnmap = NULL;
  siftop top;
  gtopology gt;
  int i, j, nri;

  switch (argc)
  {
    case 1:
      break;
    case 2:
      infn = argv[1];
      break;
    case 3:
      infn = argv[1];
      outfn = argv[2];
      break;
    default:
      fprintf(stdout, "\nError: Invalid number of generator arguments (argc=%d)\n", argc);
      return 1;
  }
  in = xfopen(infn, "rt");
  out = fopen(outfn, "wt");
  if (in == NULL)
  {
    fprintf(stderr, "\nCannot open input file %s\n\n", infn);
    return 2;
  }
  if (out == NULL)
  {
    fprintf(stderr, "\nCannot open output file %s\n\n", outfn);
    return 3;
  }

  fprintf(stdout, "\nGenerator of layered L-shape domain\n");
  fprintf(stdout, "===================================\n\n");
  
  if (cg.read(in))
  {
    xfclose(in);
    fclose(out);
    return 4;
  }

  fprintf(stdout, "\nGeneration of nodes ...");
  if (cg.gen_nodes(out))
  {
    fprintf(stdout, " FAILURE.\n");
    return 5;
  }
  else
    fprintf(stdout, " OK.\n");

  fprintf(stdout, "Generation of elements ...");
  if (cg.gen_elements(out))
  {
    fprintf(stdout, " FAILURE.\n");
    return 6;
  }
  else
    fprintf(stdout, " OK.\n");

  fprintf(stdout, "\nTopology is written in the file '%s'\n", outfn);
  fprintf(stdout, " - number of nodes in the topology: %d\n", cg.tnn);
  fprintf(stdout, " - number of elements in the topology: %d\n\n", cg.tne);

  xfclose(in);
  fclose(out);

  intf = xfopen(outfn, "rt");

  
  fprintf(stdout, "Reading of output FE topology ...");
  if (top.read(intf, 0, 1))
  {
    xfclose(in);
    fprintf(stdout, " FAILURE.\n");
    return 7;
  }
  else
    fprintf(stdout, " OK.\n");
  xfclose(intf);
  top.exporttop(&gt);

  fprintf(stdout, "Computing hanging nodes of C-profiles ...");
  cpnmap = new nodmap[cg.nn_cprof];
  for (i=cg.nn_wall, j=0; i<cg.tnn; i++, j++)
  {
    cpnmap[j].nid = i;
    cpnmap[j].x = top.nodes[i].x;
    cpnmap[j].y = top.nodes[i].y;
    cpnmap[j].z = top.nodes[i].z;
  }
  nri = cg.cprof_hang_nodes_natcoords(&gt, cg.ne_wall, cg.nn_cprof, cpnmap, 100, 1.0e-8);
  if (nri == 0)
    fprintf(stdout, " OK.\n");
  else
  {
    delete [] cpnmap;
    return 8;
  }
    
  fprintf(stdout, "Printing hanging nodes of C-profiles ...");
  out = fopen(outhnfn, "wt");
  nri = cg.print_hang_nodes(out, &gt, cg.nn_cprof, cpnmap);
  if (nri == 0)
    fprintf(stdout, " OK.\n");
  else
  {
    fprintf(stdout, " FAILURE\n");
    delete [] cpnmap;
    return 9;
  }

  fclose(out);
  delete [] cpnmap;
  return 0;
}
