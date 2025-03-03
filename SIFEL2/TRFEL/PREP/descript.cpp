#include "descript.h"
#include <stdio.h>


descript::descript()
{
  topf[0] = matf[0] = crf[0] = icf[0] = hangnf[0] = '\0';

  meshfmt = sifel;
  paral = 0;
  redgn = 0;
  matsec = crssec = no;
  matstr = crsstr = no;
  matkwd = crskwd = no;
  inicdf = no;
}



descript::~descript()
{
}



long descript::print(FILE *out)
{
  fprintf(out, "%s\n", topf);
  if (matsec == no)
    fprintf(out, "%s\n", matf);
  if (crssec == no)
    fprintf(out, "%s\n", crf);

  fprintf(out, "mesh_format %d  # mesh format indicator\n", int(meshfmt));
  fprintf(out, "edge_numbering %ld  # edge/surface property on elements indicator\n", redgn);

  if (inicdf == yes)
  {
    fprintf(out, "inicd_file %d\n", inicdf);
    fprintf(out, "%s\n", icf);
  }

  if (hangnf[0])
    fprintf(out, "hanging_nodes_file %s", hangnf);

  if ((matstr == no) || (matkwd == yes))
  {
    fprintf(out, "read_mat_strings %d\n", matstr);
    fprintf(out, "read_mat_kwd %d\n", matkwd);
  }

  if ((crsstr == no) || (crskwd == yes))
  {
    fprintf(out, "read_crs_strings %d\n", crsstr);
    fprintf(out, "read_crs_kwd %d\n", crskwd);
  }

  return 0;
}
